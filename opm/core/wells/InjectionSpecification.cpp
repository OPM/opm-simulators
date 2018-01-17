/*
  Copyright 2012  Sintef.
  Copyright 2016  Statoil.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify it under the terms
  of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  OPM is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  OPM.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdexcept>
#include "config.h"

#include <opm/common/ErrorMacros.hpp>
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
      guide_rate_(-1.0),
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
        OPM_THROW(std::domain_error, "Unknown control mode " << mode << " encountered in injection specification");
    }


    std::string
    InjectionSpecification::toString(const InjectorType& type)
    {
        switch(type) {
        case InjectorType::WATER: return "WATER";
        case InjectorType::OIL  : return "OIL"  ;
        case InjectorType::GAS  : return "GAS"  ;
        }
        OPM_THROW(std::domain_error, "Unknown injector type " << type << " encountered in injection specification");
    }


    std::string 
    InjectionSpecification::toString(const GuideRateType& type)
    {
        switch(type) {
        case GuideRateType::RAT     : return "RAT"     ;
        case GuideRateType::NONE_GRT: return "NONE_GRT";
        }
        OPM_THROW(std::domain_error, "Unknown guide rate type " << type << " encountered in injection specification");
    }
} // namespace Opm
