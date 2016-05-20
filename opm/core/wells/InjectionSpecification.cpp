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

    const std::string
    InjectionSpecification::ControlMode2String(const ControlMode& mode)
    {
        std::string stringValue;

        switch(mode) {
        case ControlMode::NONE:
            stringValue = "NONE";
            break;
        case ControlMode::RATE:
            stringValue = "RATE";
            break;
        case ControlMode::RESV:
            stringValue = "RESV";
            break;
        case ControlMode::BHP:
            stringValue = "BHP";
            break;
        case ControlMode::THP:
            stringValue = "THP";
            break;
        case ControlMode::REIN:
            stringValue = "REIN";
            break;
        case ControlMode::VREP:
            stringValue = "VREP";
            break;
        case ControlMode::GRUP:
            stringValue = "GRUP";
            break;
        case ControlMode::FLD:
            stringValue = "FLD";
            break;
        }

        return stringValue;
    }


    const std::string
    InjectionSpecification::InjectorType2String(const InjectorType& type)
    {
        std::string stringValue;

        switch(type) {
        case InjectorType::WATER:
            stringValue = "WATER";
            break;
        case InjectorType::OIL:
            stringValue = "OIL";
            break;
        case InjectorType::GAS:
            stringValue = "GAS";
            break;
        }

        return stringValue;
    }


    const std::string 
    InjectionSpecification::GuideRateType2String(const GuideRateType& type)
    {
        std::string stringValue;
        
        switch(type) {
        case GuideRateType::RAT:
            stringValue = "RAT";
            break;
        case GuideRateType::NONE_GRT:
            stringValue = "NONE_GRT";
            break;
        }

        return stringValue;
    }
} // namespace Opm
