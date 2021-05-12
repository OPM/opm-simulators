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

namespace Opm::FlowKeywordValidation
{

template<>
const KeywordValidation::PartiallySupportedKeywords<std::string>&
partiallySupported()
{
    static const KeywordValidation::PartiallySupportedKeywords<std::string> partially_supported_keywords_strings = {
        {
            "COMPORD",
            {
                {2, {false, {"INPUT"}, std::nullopt}}, // ORDER_TYPE
            },
        },
        {
            "ENDSCALE",
            {
                {1, {false, {"NODIR"}, std::nullopt}}, // DIRECT
                {2, {false, {"REVERS"}, std::nullopt}}, // IRREVERS
            },
        },
        {
            "PINCH",
            {
                {2, {false, {"GAP"}, std::nullopt}}, // GAP
                {4, {false, {"TOPBOT"}, std::nullopt}}, // PINCHOUT_OPTION
            },
        },
    };

    return partially_supported_keywords_strings;
}

template<>
const KeywordValidation::PartiallySupportedKeywords<int>&
partiallySupported()
{
    static const KeywordValidation::PartiallySupportedKeywords<int> partially_supported_keywords_int = {
        {
            "EHYSTR",
            {
                {2, {false, {0}, std::nullopt}}, //relative_perm_hyst
            },
        },
    };

    return partially_supported_keywords_int;
}

} // namespace Opm::FlowKeywordValidation
