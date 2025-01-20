/*
  Copyright 2021 Equinor.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/utils/FullySupportedFlowKeywords.hpp>

using namespace Opm::KeywordValidation;
namespace Opm::FlowKeywordValidation
{

template <>
const SupportedKeywordItems<std::string>&
fullySupported()
{
   static const SupportedKeywordItems<std::string> fully_supported_keywords_strings = {
         {
            "GEFAC",
            {
               {3,{true, is_bool_convertible {}, "GEFAC(GRPNETWK): String value must be convertible to bool."}}, // USE_GEFAC_IN_NETWORK
            },
         },
         {
            "NEXTSTEP",
            {
               {2,{true, is_bool_convertible {}, "NEXTSTEP(NSTEP2): String value must be convertible to bool."}}, // APPLY_TO_FUTURE_REPORT_STEPS
            },
         },
         {
            "WCONHIST",
            {
               {3,{true, allow_values<std::string> {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP"}, "WCONHIST(TARGET): should be set to ORAT/WRAT/GRAT/LRAT/RESV or BHP"}}, // CMODE
            },
         },
         {
            "WEFAC",
            {
               {3,{true, is_bool_convertible {}, "WEFAC(WELNETWK): String value must be convertible to bool."}}, // USE_WEFAC_IN_NETWORK
            },
         },
   };

   return fully_supported_keywords_strings;
}



template <>
const SupportedKeywordItems<int>&
fullySupported()
{
   static const SupportedKeywordItems<int> fully_supported_keywords_int = {
   };

   return fully_supported_keywords_int;
}

template <>
const SupportedKeywordItems<double>&
fullySupported()
{
   static const SupportedKeywordItems<double> fully_supported_keywords_double = {
         {
            "NEFAC",
            {
               {2,{true, [](double x) { return x > 0 && x <= 1.0; }, "NEFAC(EFF_FACTOR: Efficiency must be in the range (0,1]"}}, // NETWORK_EFF_FACTOR
            },
         },
         {
            "WPIMULT",
            {
               {2,{true, [](double x) { return x > 0; }, "WPIMULT(PIMULT): Well PI multiplier must be a positive number."}}, // PI_MULTIPLIER
            },
         },
   };

   return fully_supported_keywords_double;
}

} // namespace Opm::FlowKeywordValidation
