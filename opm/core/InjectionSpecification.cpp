#include <opm/core/InjectionSpecification.hpp>
namespace Opm
{

    InjectionSpecification::InjectionSpecification()
    : injector_type_(WATER),
      control_mode_(NONE),
      surface_flow_max_rate_(1e100),
      reinjection_fraction_target_(0.0),
      BHP_limit_(1e100),
      fluid_volume_max_rate_(1e100)
    {

    }


} // namespace Opm
