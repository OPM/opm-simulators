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

    std::string
    InjectionSpecification::toString(const ControlMode& mode)
    {
        switch(mode) {
        case ControlMode::NONE: return "NONE";
        case ControlMode::RATE: return "RATE";
        case ControlMode::RESV: return "RESV";
        case ControlMode::BHP : return "BHP" ;
        case ControlMode::THP : return "THP" ;
        case ControlMode::REIN: return "REIN";
        case ControlMode::VREP: return "VREP";
        case ControlMode::GRUP: return "GRUP";
        case ControlMode::FLD : return "FLD" ;
        }
    }


    std::string
    InjectionSpecification::toString(const InjectorType& type)
    {
        switch(type) {
        case InjectorType::WATER: return "WATER";
        case InjectorType::OIL  : return "OIL"  ;
        case InjectorType::GAS  : return "GAS"  ;
        }
    }


    std::string 
    InjectionSpecification::toString(const GuideRateType& type)
    {
        switch(type) {
        case GuideRateType::RAT     : return "RAT"     ;
        case GuideRateType::NONE_GRT: return "NONE_GRT";
        }
    }
} // namespace Opm
