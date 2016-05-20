#include "config.h"
#include <opm/core/wells/ProductionSpecification.hpp>

namespace Opm
{

    ProductionSpecification::ProductionSpecification()
    :
      control_mode_(NONE),
      procedure_(NONE_P),
      oil_max_rate_(-1e100),
      water_max_rate_(-1e100),
      gas_max_rate_(-1e100),
      liquid_max_rate_(-1e100),
      reservoir_flow_max_rate_(-1e100),
      BHP_limit_(-1e100),
      guide_rate_(1.0),
      guide_rate_type_(NONE_GRT)
    {
    }


    const std::string
    ProductionSpecification::ControlMode2String(const ControlMode& mode)
    {
        std::string stringValue;

        switch(mode) {
        case ControlMode::NONE:
            stringValue = "NONE";
            break;
        case ControlMode::ORAT:
            stringValue = "ORAT";
            break;
        case ControlMode::WRAT:
            stringValue = "WRAT";
            break;
        case ControlMode::GRAT:
            stringValue = "GRAT";
            break;
        case ControlMode::LRAT:
            stringValue = "LRAT";
            break;
        case ControlMode::CRAT:
            stringValue = "CRAT";
            break;
        case ControlMode::RESV:
            stringValue = "RESV";
            break;
        case ControlMode::PRBL:
            stringValue = "RPBL";
            break;
        case ControlMode::BHP:
            stringValue = "BHP";
            break;
        case ControlMode::THP:
            stringValue = "THP";
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
    ProductionSpecification::Procedure2String(const Procedure& type)
    {
        std::string stringValue;

        switch(type) {
        case Procedure::NONE_P:
            stringValue = "NONE_P";
            break;
        case Procedure::RATE:
            stringValue = "RATE";
            break;
        case Procedure::WELL:
            stringValue = "WELL";
            break;
        }

        return stringValue;
    }


    const std::string 
    ProductionSpecification::GuideRateType2String(const GuideRateType& type)
    {
        std::string stringValue;
        
        switch(type) {
        case GuideRateType::OIL:
            stringValue = "OIL";
            break;
        case GuideRateType::GAS:
            stringValue = "GAS";
            break;
        case GuideRateType::WATER:
            stringValue = "WATER";
            break;
        case GuideRateType::NONE_GRT:
            stringValue = "NONE_GRT";
            break;
        }

        return stringValue;
    }


}
