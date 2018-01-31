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
      guide_rate_(-1.0),
      guide_rate_type_(NONE_GRT)
    {
    }


    std::string
    ProductionSpecification::toString(const ControlMode& mode)
    {
        switch(mode) {
        case ControlMode::NONE: return "NONE";
        case ControlMode::ORAT: return "ORAT";
        case ControlMode::WRAT: return "WRAT";
        case ControlMode::GRAT: return "GRAT";
        case ControlMode::LRAT: return "LRAT";
        case ControlMode::CRAT: return "CRAT";
        case ControlMode::RESV: return "RESV";
        case ControlMode::PRBL: return "RPBL";
        case ControlMode::BHP : return "BHP" ;
        case ControlMode::THP : return "THP" ;
        case ControlMode::GRUP: return "GRUP";
        case ControlMode::FLD : return "FLD" ;
        }
        OPM_THROW(std::domain_error, "Unknown control mode " << mode << " encountered in production specification");
    }


    std::string 
    ProductionSpecification::toString(const Procedure& type)
    {
        switch(type) {
        case Procedure::NONE_P: return "NONE_P";
        case Procedure::RATE  : return "RATE"  ;
        case Procedure::WELL  : return "WELL"  ;
        }
        OPM_THROW(std::domain_error, "Unknown procedure " << type << " encountered in production specification");
    }


    std::string 
    ProductionSpecification::toString(const GuideRateType& type)
    {
        switch(type) {
        case GuideRateType::OIL     : return "OIL"     ;
        case GuideRateType::GAS     : return "GAS"     ;
        case GuideRateType::WATER   : return "WATER"   ;
        case GuideRateType::LIQ     : return "LIQ"     ;
        case GuideRateType::NONE_GRT: return "NONE_GRT";
        }
        OPM_THROW(std::domain_error, "Unknown guide rate type " << type << " encountered in production specification");
    }


}
