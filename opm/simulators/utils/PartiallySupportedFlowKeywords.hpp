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

#ifndef OPM_PARTIALLYSUPPORTEDFLOWKEYWORDS_HEADER_INCLUDED
#define OPM_PARTIALLYSUPPORTEDFLOWKEYWORDS_HEADER_INCLUDED


#include <map>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/simulators/flow/KeywordValidation.hpp>

/*
    Here keywords are defined that are supported by flow, but have items that
    are only partially supported.

    The keywords are specified in a mapping with the keyword names as keys, and
    values that describe the set of supported items. These are described by a
    mapping from the item name to a struct of properties, defined in KeywordValidation.hpp.

    This struct has the following fields:

    critical (bool) : if this is a critical error.
    permitted_values (vector of strings) : the list of values that is allowed.
    message (itemal string): an optional message to add to the error reported by flow.

    Below is the set of partiall supported keywords, currently used by flow.
*/

namespace Opm::FlowKeywordValidation
{


const Opm::KeywordValidation::PartiallySupportedKeywords<std::string> partially_supported_keywords_strings = {
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

const Opm::KeywordValidation::PartiallySupportedKeywords<int> partially_supported_keywords_int = {
    {
        "EHYSTR",
        {
            {2, {false, {0}, std::nullopt}}, //relative_perm_hyst
        },
    },
};


} // namespace Opm::FlowKeywordValidation


#endif
