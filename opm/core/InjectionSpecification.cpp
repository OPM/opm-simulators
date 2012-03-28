#include <opm/core/InjectionSpecification.hpp>
namespace Opm {

InjectionSpecification::InjectionSpecification()
: component_(WATER), control_mode_(NONE), surface_injection_target_(0.0),
reinjection_fraction_target_(0.0), BHP_target_(0.0) {

}

InjectionSpecification::~InjectionSpecification() {
}

InjectionSpecification::Component InjectionSpecification::component() {
    return component_;
}

void InjectionSpecification::set_component(InjectionSpecification::Component comp) {
    component_ = comp;
}

InjectionSpecification::ControlMode InjectionSpecification::control_mode() {
    return control_mode_;
}

void InjectionSpecification::set_control_mode(InjectionSpecification::ControlMode mode) {
    control_mode_ = mode;
}

/// \returns 0 if no limit, else the target/limit of the surface Injector 
///          rate.
double InjectionSpecification::surface_injection_target() {
    return surface_injection_target_;
}

void InjectionSpecification::set_surface_injection_target(double target) {
    surface_injection_target_ = target;
}

/// \returns 0 if no limit, else the target/limit of the reInjector fraction

double InjectionSpecification::reinjection_fraction_target() {
    return reinjection_fraction_target_;
}

void InjectionSpecification::set_reinjection_fraction_target(double target) {
    reinjection_fraction_target_ = target;
}

/// \returns 0 if no limit, else the target/limit of the BHP

double InjectionSpecification::BHP_target() {
    return BHP_target_;
}

void InjectionSpecification::set_BHP_target(double target) {
    BHP_target_ = target;
}

} // namespace Opm