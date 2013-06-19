#include "config.h"
#include <opm/core/wells/InjectionSpecification.hpp>

namespace Opm
{

    InjectionSpecification::InjectionSpecification()
    : injector_type_(WATER),
      control_mode_(NONE),
      surface_flow_max_rate_(-1e100),
      reservoir_flow_max_rate_(-1e100),
      BHP_limit_(-1e100),
      reinjection_fraction_target_(1),
      voidage_replacment_fraction_(1),
      guide_rate_(1.0),
      guide_rate_type_(NONE_GRT)
    {

    }

} // namespace Opm
