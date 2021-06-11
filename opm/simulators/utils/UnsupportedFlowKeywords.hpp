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

#ifndef OPM_UNSUPPORTEDFLOWKEYWORDS_HEADER_INCLUDED
#define OPM_UNSUPPORTEDFLOWKEYWORDS_HEADER_INCLUDED


#include <opm/simulators/flow/KeywordValidation.hpp>

/*
    Here the keywords are defined that are not suppored by flow. When parsing
    the deck, unsupported keywords will be flagged and reported. It is possible
    to mark a keyword as 'critical'. In this case an error will be raised and
    flow will not proceed when it encounters this keyword.

    The unsupported keywords are specified by a mapping from the keyword name to
    a struct with properties of that keyword. The struct is defined in
    KeywordValidation.hpp and contains the following fields:

    critical (bool) : set to true if the keywords is considered critical.
    message (optional string): set to an optional message string that is added
                               to the error.

    Below is the std::map with the keywords that are not supported by flow.

    Most entries in the map below are not critical ('false') and do not have an
    additional message. Set the first entry to 'true' to make a keyword
    critical. A message can be added by replacing std::nullopt by a string
    literal.
*/


namespace Opm::FlowKeywordValidation
{

const KeywordValidation::UnsupportedKeywords& unsupportedKeywords();

} // namespace Opm::FlowKeywordValidation


#endif
