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
const SupportedKeywords<std::string>&
fullySupported()
{
   static const SupportedKeywords<std::string> fully_supported_keywords_strings = {
         {
            "WEFAC",
            {
               {3,{true, is_bool_convertible {}, "WEFAC(WELNETWK): String value must be convertible to bool."}}, // USE_WEFAC_IN_NETWORK
            },
         },
         {
            "GEFAC",
            {
               {3,{true, is_bool_convertible {}, "GEFAC(GRPNETWK): String value must be convertible to bool."}}, // USE_GEFAC_IN_NETWORK
            },            
         },
   };

   return fully_supported_keywords_strings;
}



template <>
const SupportedKeywords<int>&
fullySupported()
{
   static const SupportedKeywords<int> fully_supported_keywords_int = {
   };

   return fully_supported_keywords_int;
}

template <>
const SupportedKeywords<double>&
fullySupported()
{
   static const SupportedKeywords<double> fully_supported_keywords_double = {
         {
            "WEFAC",
            {
               {2,{true, [](double x) { return x > 0; }, "WEFAC(FACTOR): Well efficiency must be a positive number."}}, // WELL_EFFICIENCY
            },
         },
         {
            "GEFAC",
            {
               {2,{true, [](double x) { return x > 0; }, "GEFAC(FACTOR): Well efficiency must be a positive number."}}, // WELL_EFFICIENCY
            },            
         },
   };

   return fully_supported_keywords_double;
}

} // namespace Opm::FlowKeywordValidation
