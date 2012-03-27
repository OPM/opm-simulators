#include <opm/core/ProductionSpecification.hpp>

namespace Opm { 
ProductionSpecification::ProductionSpecification() : 
        component_(OIL),
        control_mode_(NONE_CM),
        procedure_(NONE_P),
        oil_production_target_(-1.0),
        water_production_target_(-1.0),
        liquid_production_target_(-1.0),
        BHP_target_(-1.0)
{
}

ProductionSpecification::~ProductionSpecification() {
}

void ProductionSpecification::set_BHP_target(double BHP_target_) {
    this->BHP_target_ = BHP_target_;
}

double ProductionSpecification::get_BHP_target() const {
    return BHP_target_;
}

void ProductionSpecification::set_liquid_production_target(double liquid_production_target_) {
    this->liquid_production_target_ = liquid_production_target_;
}

double ProductionSpecification::get_liquid_production_target() const {
    return liquid_production_target_;
}

void ProductionSpecification::set_water_production_target(double water_production_target_) {
    this->water_production_target_ = water_production_target_;
}

double ProductionSpecification::get_water_production_target() const {
    return water_production_target_;
}

void ProductionSpecification::set_oil_production_target(double oil_production_target_) {
    this->oil_production_target_ = oil_production_target_;
}

double ProductionSpecification::get_oil_production_target() const {
    return oil_production_target_;
}

void ProductionSpecification::set_procedure(ProductionSpecification::Procedure procedure_) {
    this->procedure_ = procedure_;
}

ProductionSpecification::Procedure ProductionSpecification::get_procedure() const {
    return procedure_;
}

void ProductionSpecification::set_control_mode(ProductionSpecification::ControlMode control_mode_) {
    this->control_mode_ = control_mode_;
}

ProductionSpecification::ControlMode ProductionSpecification::get_control_mode() const {
    return control_mode_;
}

void ProductionSpecification::set_component(ProductionSpecification::Component component_) {
    this->component_ = component_;
}

ProductionSpecification::Component ProductionSpecification::get_component() const {
    return component_;
}
}
