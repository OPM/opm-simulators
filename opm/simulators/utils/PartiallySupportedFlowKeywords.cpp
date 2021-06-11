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

#include <opm/simulators/utils/PartiallySupportedFlowKeywords.hpp>

using namespace Opm::KeywordValidation;
namespace Opm::FlowKeywordValidation
{

template <>
const PartiallySupportedKeywords<std::string>&
partiallySupported()
{
    static const PartiallySupportedKeywords<std::string> partially_supported_keywords_strings = {
        {
            "COMPORD",
            {
                {2, {false, allow_values<std::string>{"INPUT"}, std::nullopt}}, // ORDER_TYPE
            },
        },
        {
            "ENDSCALE",
            {
                {1, {false, allow_values<std::string>{"NODIR"}, std::nullopt}}, // DIRECT
                {2, {false, allow_values<std::string>{"REVERS"}, std::nullopt}}, // IRREVERS
            },
        },
        {
            "PINCH",
            {
                {2, {false, allow_values<std::string>{"GAP"}, std::nullopt}}, // GAP
                {4, {false, allow_values<std::string>{"TOPBOT"}, std::nullopt}}, // PINCHOUT_OPTION
            },
        },
    };

    return partially_supported_keywords_strings;
}

template <>
const PartiallySupportedKeywords<int>&
partiallySupported()
{
    static const PartiallySupportedKeywords<int> partially_supported_keywords_int = {
        {
            "EHYSTR",
            {
                {2, {false, allow_values<int>{0}, std::nullopt}}, // relative_perm_hyst
            },
        },
    };

    return partially_supported_keywords_int;
}

template <>
const PartiallySupportedKeywords<double>&
partiallySupported()
{
    static const PartiallySupportedKeywords<double> partially_supported_keywords_double = {};

    return partially_supported_keywords_double;
}

} // namespace Opm::FlowKeywordValidation
