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


#include <opm/simulators/flow/KeywordValidation.hpp>

#include <string>

/*
    Here keywords are defined that are supported by flow, but have items that
    are only partially supported.

    The keywords are specified in a mapping with the keyword names as keys, and
    values that describe the set of supported items. These are described by a
    mapping from the item name to a struct of properties, defined in KeywordValidation.hpp.

    This struct has the following fields:

    critical (bool) : if this is a critical error.
    validator (function wrapper) : A function wrapper object that is used to test values.
    message (itemal string): an optional message to add to the error reported by flow.

    For convenience there is a small class KeywordValidation::allow_values which
    can be initialized with a list of permitted values, and used as a validator.
*/

namespace Opm::FlowKeywordValidation
{

template <typename T>
const KeywordValidation::PartiallySupportedKeywords<T>& partiallySupported();

} // namespace Opm::FlowKeywordValidation


#endif
