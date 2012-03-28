#include <opm/core/ProductionSpecification.hpp>

namespace Opm
{

    ProductionSpecification::ProductionSpecification()
    : component_(OIL),
      control_mode_(NONE_CM),
      procedure_(NONE_P),
      oil_production_target_(-1.0),
      water_production_target_(-1.0),
      liquid_production_target_(-1.0),
      BHP_target_(-1.0)
    {
    }

}
