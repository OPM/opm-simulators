#include <opm/core/InjectionSpecification.hpp>
namespace Opm
{

    InjectionSpecification::InjectionSpecification()
    : component_(WATER), control_mode_(NONE), surface_injection_target_(0.0),
    reinjection_fraction_target_(0.0), BHP_target_(0.0)
    {

    }


} // namespace Opm