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

#ifndef VALIDATION_FUNCTIONS_HPP
#define VALIDATION_FUNCTIONS_HPP

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm {
class DeckKeyword;
}

namespace Opm::KeywordValidation {

struct ValidationError;

using ValidationFunction = std::function<void(const DeckKeyword&,
                                              std::vector<ValidationError>&)>;

// This is a mapping between keyword names and small functions
// for validation of special keywords.
std::unordered_map<std::string, ValidationFunction>
specialValidation();

} // namespace Opm::KeywordValidation

#endif // VALIDATION_FUNCTIONS_HPP
